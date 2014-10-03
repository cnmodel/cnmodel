"""
Utility class for spawning and controlling CLI processes. 
See: http://stackoverflow.com/questions/375427/non-blocking-read-on-a-subprocess-pipe-in-python
"""
import sys
from subprocess import PIPE, Popen
from threading  import Thread

try:
    from Queue import Queue, Empty
except ImportError:
    from queue import Queue, Empty  # python 3.x

ON_POSIX = 'posix' in sys.builtin_module_names


class Process(object):
    """
    Process encapsulates a subprocess with queued stderr/stdout pipes. 
    Wraps most methods from subprocess.Popen.
    
    For non-blocking reads, use proc.stdout.get_nowait() or 
    proc.stdout.get(timeout=0.1).
    """
    def __init__(self, exec_args):
        self.proc = Popen(exec_args, stdout=PIPE, stdin=PIPE, stderr=PIPE, 
                          bufsize=1, close_fds=ON_POSIX, universal_newlines=True)
        self.stdin = self.proc.stdin
        #self.stdout = self.proc.stdout
        #self.stderr = self.proc.stderr
        self.stdout = PipeQueue(self.proc.stdout)
        self.stderr = PipeQueue(self.proc.stderr)
        for method in ['poll', 'wait', 'send_signal', 'kill', 'terminate']:
            setattr(self, method, getattr(self.proc, method))


class PipeQueue(Queue):
    """
    Queue that starts a second process to monitor a PIPE for new data.
    """
    def __init__(self, pipe):
        Queue.__init__(self)
        self.thread = Thread(target=self.enqueue_output, args=(pipe, self))
        self.thread.daemon = True # thread dies with the program
        self.thread.start()
        
    @staticmethod
    def enqueue_output(out, queue):
        for line in iter(out.readline, b''):
            queue.put(line)
            #print "READ: " + repr(line)
        out.close()

    def read(self):
        """
        Read all available lines from the queue, concatenated into a single string.
        """
        out = ""
        while True:
            try:
                out += self.get_nowait()
            except Empty:
                break
        return out
        
    def readline(self, timeout=None):
        return self.get(timeout=timeout)

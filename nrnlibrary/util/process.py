"""
Utility class for spawning and controlling CLI processes. 
See: http://stackoverflow.com/questions/375427/non-blocking-read-on-a-subprocess-pipe-in-python
"""
import sys, time
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
    def __init__(self, exec_args, cwd=None):
        self.proc = Popen(exec_args, stdout=PIPE, stdin=PIPE, stderr=PIPE, 
                          bufsize=1, close_fds=ON_POSIX, universal_newlines=True,
                          cwd=cwd)
        self.stdin = self.proc.stdin
        #self.stdin = PipePrinter(self.proc.stdin)
        self.stdout = PipeQueue(self.proc.stdout)
        self.stderr = PipeQueue(self.proc.stderr)
        for method in ['poll', 'wait', 'send_signal', 'kill', 'terminate']:
            setattr(self, method, getattr(self.proc, method))


class PipePrinter(object):
    """ For debugging writes to a pipe.
    """
    def __init__(self, pipe):
        self._pipe = pipe
        
    def __getattr__(self, attr):
        return getattr(self._pipe, attr)
    
    def write(self, strn):
        print "WRITE:" + repr(strn)
        return self._pipe.write(strn)


class PipeQueue(Queue):
    """
    Queue that starts a second process to monitor a PIPE for new data.
    This is needed to allow non-blocking pipe reads.
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
        """ Read a single line from the queue.
        """
        # we break this up into multiple short reads to allow keyboard 
        # interrupts
        start = time.time()
        ret = ""
        while True:
            if timeout is not None:
                remaining = start + timeout - time.time()
                if remaining <= 0:
                    return ""
            else:
                remaining = 1
            
            try:
                return self.get(timeout=min(0.1, remaining))
            except Empty:
                pass

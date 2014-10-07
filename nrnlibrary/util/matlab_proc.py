"""
Simple system for interfacing with a MATLAB process using stdin/stdout pipes.

"""
from process import Process
from StringIO import StringIO
import scipy.io
import numpy as np
import tempfile
import os

class MatlabProcess(Process):
    """ This class starts a new matlab process, allowing remote control of
    the interpreter and transfer of data between python and matlab.
    """
    
    # Implements a MATLAB CLI that is a bit easier for us to parse
    _bootstrap = r"""
    while true
        fprintf('\n::ready\n');
        
        % Accumulate command lines until we see ::cmd_done
        cmd = '';
        empty_count = 0;
        while true
            line = input('', 's');
            if length(line) == 0
                % If we encounter many empty lines, assume the host process
                % has ended. Is there a better way to detect this?
                empty_count = empty_count + 1;
                if empty_count > 100
                    fprintf('Shutting down!\n');
                    exit;
                end
            end
            if strcmp(line, '::cmd_done')
                break
            end
            cmd = [cmd, line, sprintf('\n')]
        end
        
        % Evaluate command
        try
            %fprintf('EVAL: %s\n', cmd);
            eval(cmd);
            fprintf('\n::ok\n');
        catch err
            fprintf('\n::err\n');
            fprintf(['::message:', err.message, '\n']);
            fprintf(['::identifier:', err.identifier, '\n']);
            for i = 1:length(err.stack)
                frame = err.stack(i,1);
                fprintf(['::stack:', frame.name, ' in ', frame.file, ' line ', frame.line]);
            end
        end
    end
    """
    
    def __init__(self, executable='matlab', **kwds):
        self.cmd_index = 0
        Process.__init__(self, [executable, '-nodesktop', '-nosplash'], **kwds)
        # Wait a moment for MATLAB to start up, 
        # read the version string
        while True:
            line = self.stdout.readline()
            if 'Copyright' in line:
                # next line is version info
                self.version_str = self.stdout.readline().strip()
                break
            
        # start input loop
        self.stdin.write(self._bootstrap)
        
        # wait for input loop to be ready
        while True:
            line = self.stdout.readline()
            if line == '::ready\n':
                break
            
        
    def __call__(self, cmd, timeout=5.0):
        """
        Execute the specified statement(s) on the MATLAB interpreter and return
        the output string or raise an exception if there was an error.
        """
        if cmd[-1] != '\n':
            cmd += '\n'
        cmd += "::cmd_done\n"
        self.stdout.read()
        self.stdin.write(cmd)
        
        return self._parse_result()
    
    def _parse_result(self):
        output = []
        while True:
            line = self.stdout.readline()
            if line == '::ready\n':
                break
            output.append(line)
                
        for i in reversed(range(len(output))):
            line = output[i]
            if line == '::ok\n':
                return ''.join(output[:i])
            elif line == '::err\n':
                raise MatlabError(output[i+1:])
            
        raise RuntimeError("No success/failure code found in output (printed above).")

    def get(self, name):
        """
        Transfer an object from MATLAB to Python.
        """
        assert isinstance(name, str)
        tmp = tempfile.mktemp(suffix='.mat')
        out = self("save('%s', '%s', '-v7')" % (tmp, name))
        objs = scipy.io.loadmat(tmp)
        os.remove(tmp)
        return objs[name]

    def set(self, **kwds):
        """
        Transfer an object from Python to MATLAB and assign it to the given
        variable name.
        """
        tmp = tempfile.mktemp(suffix='.mat')
        scipy.io.savemat(tmp, kwds)
        self("load('%s')" % tmp)
        os.remove(tmp)
                
    def get_via_pipe(self, name):
        """
        Transfer an object from MATLAB to Python.
        
        This method sends data over the pipe, but is less reliable than get().
        """
        assert isinstance(name, str)
        out = self("save('stdio', '%s', '-v7')" % name)
        start = stop = None
        for i, line in enumerate(out):
            if line.startswith('start_binary'):
                start = i
            elif line.startswith('::ok'):
                stop = i
        data = ''.join(out[start+1:stop])
        io = StringIO(data[:-1])
        objs = scipy.io.loadmat(io)
        return objs[name]

    def set_via_pipe(self, **kwds):
        """
        Transfer an object from Python to MATLAB and assign it to the given
        variable name.
        
        This method sends data over the pipe, but is less reliable than set().
        """
        io = StringIO()
        scipy.io.savemat(io, kwds)
        io.seek(0)
        strn = io.read()
        self.stdout.read()
        self.stdin.write("load('stdio')\n::cmd_done\n")
        while True:
            line = self.stdout.readline()
            if line == 'ack load stdio\n':
                # now it is safe to send data
                break
        self.stdin.write(strn)
        self.stdin.write('\n')
        while True:
            line = self.stdout.readline()
            if line == 'ack load finished\n':
                break
        self._parse_result()

    def __getattr__(self, name):
        return MatlabCallable(self, name)
        

class MatlabCallable(object):
    """
    Proxy to a MATLAB function. 
    
    Calling this object transfers the arguments to MATLAB, invokes the function
    remotely, and then transfers the return value back to Python.
    """
    def __init__(self, proc, name):
        self._proc = proc
        self._name = name
        self._nargout = None

    @property
    def nargout(self):
        """ Number of output arguments for this function.
        
        For some functions, requesting nargout() will fail. In these cases,
        the nargout property must be set manually before calling the function.
        """
        if self._nargout is None:
            nargvar = "%s_nargs_%d" % (self._name, id(self))
            try:
                self._proc("%s = nargout('%s');" % (nargvar, self._name))
                self._nargout = self._proc.get(nargvar)
            finally:
                self._proc("clear %s" % nargvar)
        return self._nargout
    
    @nargout.setter
    def nargout(self, n):
        self._nargout = n
        
    def __call__(self, *args):
        # store args to temporary variables
        argnames = ['%s_%d_%d' % (self._name, i, id(self)) for i in range(len(args))]
        args = dict(zip(argnames, args))
        self._proc.set(**args)
        
        try:
            # get number of output args
            nargs = self.nargout
        
            # invoke function, fetch return value(s)
            retvars = ["%s_rval_%d_%d" % (self._name, i, id(self)) for i in range(nargs)]
            cmd = "[%s] = %s(%s)" % (','.join(retvars), self._name, ','.join(argnames))
            self._proc(cmd)
            ret = [self._proc.get(var) for var in retvars]
            if len(ret) == 1:
                ret = ret[0]
            return ret
        
        finally:
            # clear all temp variables
            cmd = "clear %s" % (' '.join(argnames + retvars))
            self._proc(cmd)
        
        return ret
        

class MatlabError(Exception):
    def __init__(self, output):
        for line in output:
            self.stack = []
            if line.startswith('::message:'):
                self.message = line[10:].strip()
            elif line.startswith('::identifier:'):
                self.identifier = line[13:].strip()
            elif line.startswith('::stack:'):
                self.stack.append(line[8:].strip())

    def __repr__(self):
        return "MatlabError(message=%s, identifier=%s)" % (repr(self.message), repr(self.identifier))
    
    def __str__(self):
        if len(self.stack) > 0:
            stack = "\nMATLAB Stack:\n%s\nMATLAB Error: " % '\n'.join(self.stack) 
            return stack + self.message
        else:
            return self.message


if __name__ == '__main__':
    p = MatlabProcess()
    io = StringIO()
    scipy.io.savemat(io, {'x': 1})
    io.seek(0)
    strn = io.read()
    

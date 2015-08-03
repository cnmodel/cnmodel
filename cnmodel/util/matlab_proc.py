"""
Simple system for interfacing with a MATLAB process using stdin/stdout pipes.

"""
from process import Process
from StringIO import StringIO
import scipy.io
import numpy as np
import tempfile
import os, sys, glob, signal


class MatlabProcess(object):
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
            cmd = [cmd, line, sprintf('\n')];
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
    
    def __init__(self, executable=None, **kwds):
        self.__refs = {}
        
        # Decide which executables to try
        if executable is not None:
            execs = [executable]
        else:
            execs = ['matlab']  # always pick the matlab in the path, if available
            if sys.platform == 'darwin':
                installed = glob.glob('/Applications/MATLAB_R*')
                installed.sort(reverse=True)
                execs.extend([os.path.join(p, 'bin', 'matlab') for p in installed])
                
        # try starting each in order until one works
        self.__proc = None
        for exe in execs:
            try:
                self.__proc = Process([exe, '-nodesktop', '-nosplash'], **kwds)
                break
            except Exception as e:
                last_exception = e
                pass
            
        # bail out if we couldn't start any
        if self.__proc is None:
            raise RuntimeError("Could not start MATLAB process.\nPaths attempted: %s "
                               "\nLast error: %s" % (str(execs), str(e)))

        
        # Wait a moment for MATLAB to start up, 
        # read the version string
        while True:
            line = self.__proc.stdout.readline()
            if 'Copyright' in line:
                # next line is version info
                self.__version_str = self.__proc.stdout.readline().strip()
                break
            
        # start input loop
        self.__proc.stdin.write(self._bootstrap)
        
        # wait for input loop to be ready
        while True:
            line = self.__proc.stdout.readline()
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
        self.__proc.stdout.read()
        self.__proc.stdin.write(cmd)
        
        return self._parse_result()
    
    def _parse_result(self):
        output = []
        while True:
            line = self.__proc.stdout.readline()
            if line == '::ready\n':
                break
            output.append(line)
                
        for i in reversed(range(len(output))):
            line = output[i]
            if line == '::ok\n':
                return ''.join(output[:i])
            elif line == '::err\n':
                raise MatlabError(output[i+1:], output[:i])
            
        raise RuntimeError("No success/failure code found in output (printed above).")

    def _get(self, name):
        """
        Transfer an object from MATLAB to Python.
        """
        assert isinstance(name, str)
        tmp = tempfile.mktemp(suffix='.mat')
        out = self("save('%s', '%s', '-v7')" % (tmp, name))
        objs = scipy.io.loadmat(tmp)
        os.remove(tmp)
        return objs[name]

    def _set(self, **kwds):
        """
        Transfer an object from Python to MATLAB and assign it to the given
        variable name.
        """
        tmp = tempfile.mktemp(suffix='.mat')
        scipy.io.savemat(tmp, kwds)
        self("load('%s')" % tmp)
        os.remove(tmp)
                
    def _get_via_pipe(self, name):
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

    def _set_via_pipe(self, **kwds):
        """
        Transfer an object from Python to MATLAB and assign it to the given
        variable name.
        
        This method sends data over the pipe, but is less reliable than set().
        """
        io = StringIO()
        scipy.io.savemat(io, kwds)
        io.seek(0)
        strn = io.read()
        self.__proc.stdout.read()
        self.__proc.stdin.write("load('stdio')\n::cmd_done\n")
        while True:
            line = self.__proc.stdout.readline()
            if line == 'ack load stdio\n':
                # now it is safe to send data
                break
        self.__proc.stdin.write(strn)
        self.__proc.stdin.write('\n')
        while True:
            line = self.__proc.stdout.readline()
            if line == 'ack load finished\n':
                break
        self._parse_result()

    def exist(self, name):
        if name == 'exist':
            return 5
        else:
            for line in self('exist %s' % name).split('\n'):
                try:
                    return int(line.strip())
                except ValueError:
                    pass

    def __getattr__(self, name):
        if name not in self.__refs:
            ex = self.exist(name)
            if ex == 0:
                raise AttributeError("No object named '%s' in matlab workspace." % name)
            if ex in (2, 3, 5):
                r = MatlabFunction(self, name)
            elif ex == 1:
                r = self._mkref(name)
            self.__refs[name] = r
        return self.__refs[name]
        
    def _mkref(self, name):
        assert name not in self.__refs
        ref = MatlabReference(self, name)
        self.__refs[name] = ref
        return ref
    
    def __setattr__(self, name, value):
        if name.startswith('_MatlabProcess__'):
            object.__setattr__(self, name, value)
        else:
            self._set(**{name:value})


class MatlabReference(object):
    """ Reference to a variable in the matlab workspace.
    """
    def __init__(self, proc, name):
        self._proc = proc
        self._name = name

    @property
    def name(self):
        return self._name

    def get(self):
        return self._proc._get(self._name)

    def clear(self):
        self._proc("clear %s;" % self._name)


class MatlabFunction(object):
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
            cmd = "fprintf('%%d\\n', nargout('%s'));" % (self._name)
            ret = self._proc(cmd)
            self._nargout = int(ret.strip())
        return self._nargout
    
    @nargout.setter
    def nargout(self, n):
        self._nargout = n
        
    def __call__(self, *args, **kwds):
        """
        Call this function with the given arguments.
        
        If _transfer is False, then the return values are left in MATLAB and 
        references to these values are returned instead.
        """
        import pyqtgraph as pg
        
        _transfer = kwds.pop('_transfer', True)
        assert len(kwds) == 0
        
        # for generating unique variable names
        rand = np.random.randint(1e12)
        
        # store args to temporary variables, excluding those already present
        # in the workspace
        argnames = []
        upload = {}
        for i, arg in enumerate(args):
            if isinstance(arg, MatlabReference):
                argnames.append(arg.name)
            elif np.isscalar(arg):
                argnames.append(repr(arg))
            else:
                argname = '%s_%d_%d' % (self._name, i, rand)
                argnames.append(argname)
                upload[argname] = arg
        if len(upload) > 0:
            self._proc._set(**upload)
        
        try:
            # get number of output args
            nargs = self.nargout
        
            # invoke function, fetch return value(s)
            retvars = ["%s_rval_%d_%d" % (self._name, i, rand) for i in range(nargs)]
            cmd = "[%s] = %s(%s);" % (','.join(retvars), self._name, ','.join(argnames))
            self._proc(cmd)
            if _transfer:
                ret = [self._proc._get(var) for var in retvars]
            else:
                ret = [self._proc._mkref(name) for name in retvars]
            if len(ret) == 1:
                ret = ret[0]
            else:
                ret = tuple(ret)
            return ret
        finally:
            # clear all temp variables
            clear = list(upload.keys())
            #if _transfer:
                #clear += retvars
            if len(clear) > 0:
                cmd = "clear %s;" % (' '.join(clear))
                self._proc(cmd)
        return ret
        

class MatlabError(Exception):
    def __init__(self, error, output):
        self.output = ''.join(output)
        for line in error:
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
    

# -*- encoding: utf-8 -*-
from __future__ import print_function
from collections import OrderedDict
import re


# Unified collection point for all empirically-determined biophysical
# values. Each value is a tuple (val, source). 
DATA = OrderedDict()

def get(*args, **kwds):
    """ Get a single value from the database using the supplied arguments
    to query. 
    
    Optionally, one keyword argument may be a list of values, in which case
    a dict will be returned containing {listval: dbval} pairs for each value in
    the list.
    """
    return _lookup(0, *args, **kwds)
    
def get_source(*args, **kwds):
    """ Get the source of a single value from the database using the supplied 
    arguments to query.
    
    Optionally, one keyword argument may be a list of values, in which case
    a dict will be returned containing {listval: dbval} pairs for each value in
    the list.
    """
    return _lookup(1, *args, **kwds)

def print_table(table):
    for k in DATA.keys():
        if table == k[0]:
            print('data key: ', k)
            print( DATA[k][0])

def get_table_info(table):
    """
    Return a dictionary of row and column names in the table
    """
    tinfo = {}
    for k in DATA.keys():
        if table == k[0]:
            for p in k:
                if not isinstance(p, tuple):
                    continue
                if p[0] not in tinfo.keys():
                    tinfo[p[0]] = []
                if p[1] not in tinfo[p[0]]:
                    tinfo[p[0]].append(p[1])
    return tinfo
            

def _lookup(ind, *args, **kwds):
    key = mk_key(*args, **kwds)
    if isinstance(key, dict):
        data = {}
        for k,key in key.items():
            data[k] = DATA[key][ind]
        return data
    else:
        return DATA[key][ind]
            
def setval(val, *args, **kwds):
    key = mk_key(*args, **kwds)
    oldval = None
    #change_flag = False
    if key in DATA:
#        change_flag = True  # any attempt to change key will set this
        oldval = DATA[key]  # save the previous stored value
    #     raise RuntimeError("Data key '%s' has already been set." % str(key))
    DATA[key] = val
    return oldval

def mk_key(*args, **kwds):
    # Make a unique key (or list of keys) used to access values from the 
    # database. The generated key is independent of the order that arguments
    # are specified.
    # 
    # Optionally, one keyword argument may have a list of values, in which case
    # the function will return a dict containing {listval: key} pairs for each
    # value in the list.
    listkey = None
    for k,v in kwds.items():
        if isinstance(v, (list, tuple)):
            if listkey is not None:
                raise TypeError("May only specify a list of values for one key.")
            listkey = k

    if listkey is None:
        return _mk_key(*args, **kwds)
    else:
        keys = {}
        for v in kwds[listkey]:
            kwds[listkey] = v
            keys[v] = _mk_key(*args, **kwds)
        return keys
        
def _mk_key(*args, **kwds):
    key = list(args) + list(kwds.items())
    key.sort(key=lambda a: a[0] if isinstance(a, tuple) else a)
    return tuple(key)


def add_table_data(name, row_key, col_key, data, **kwds):
    """
    Read data like::
    
        Description
        
        ------------------------------------
                col1      col2         col3
        row1    1.2  [1]  0.9e-6 [1]   27 [2]
        row2    1.7  [1]               [3]
        row3    0.93 [2]  0.3e-6       3 [2]
        
        ------------------------------------
        
        [1] citation 1
        [2] citation 2
        [3] missing because.
    
    
    """
    if isinstance(data, str) and '\xc2' in data:
        raise TypeError('Data table <%s> appears to contain unicode characters but'
                        'was not defined as unicode.' % name)
    
    lines = data.split('\n')
    
    # First, split into description, table, and sources using ----- lines
    desc = []
    table = []
    while lines:
        line = lines.pop(0)
        #print ">", line
        if re.match(r'\s*-+\s*$', line):
            #print "match!"
            break
        desc.append(line)
    while lines:
        line = lines.pop(0)
        #print ">", line
        if re.match(r'\s*-+\s*$', line):
            #print "match!"
            break
        table.append(line)
    
    #print desc 
    #print table
    
    # parse remaining lines as sources
    sources = parse_sources(lines)
    #print sources
    
    #
    # parse table
    # table might be empty, so take care of that first.
    if table == []:
        return []  # no changes

    while len(table[0].strip()) == 0:
            table.pop(0)
    
    spaces = [c == ' ' for c in table[0]]
    cols = [0] + [i for i in range(1, len(spaces)) if spaces[i-1] and not spaces[i]]
    cols = cols + [max(map(len, table))+1]
    #print spaces
    #print cols
    # Make sure columns are obeyed strictly
    for i,line in enumerate(table):
        for j, c in enumerate(cols[1:]):
            if len(line) < c:
                continue
            if line[c-1] != " ":
                print('line : ', line)
                raise Exception("Table <%s> line %d column %s does not obey column boundaries." % (name, (i, j)))
            
    # Break table into cells
    cells = []
    for line in table:
        if line.strip() != '':
            cells.append([line[cols[i]:cols[i+1]].strip() for i in range(len(cols)-1)])
    #print cells
    
    # Extract row/column names
    col_names = cells.pop(0)[1:]
    row_names = [cells[i].pop(0) for i in range(len(cells))]
    if len(set(row_names)) != len(row_names):
        for n in set(row_names):
            row_names.remove(n)
        raise NameError('Duplicate row names: %s' % row_names)
    
    # Parse cell values
    for i in range(len(cells)):
        for j in range(len(cells[0])):
            cell = cells[i][j].strip()
            m = re.match(r'([^\[]*)(\[([^\]]+)\])?', cell)  # match like "0.7 [3]"
            if m is None:
                raise ValueError("Table cell (%d, %d) has bad format: '%s'" % (i, j, cell))
            
            # parse value
            # If the value contains '±' then a tuple is returned containing the values
            # on either side.
            val, _, source = m.groups()
            #val = unicode(val)  # python 2
            val = str(val)  # python 3
            if val.strip() == '':
                val = None
            else:
                parts = val.split(u'±')
                vals = []
                for p in parts:
                    try:
                        p = int(p)
                    except ValueError:
                        try:
                            p = float(p)
                        except ValueError:
                            try:
                                 p = str(p.strip())  # allow strings to identify mechanisms also
                            except ValueError:
                                raise ValueError("Table cell (%d, %d) value has bad format: '%s'" % (i, j, val))
                    vals.append(p)
                if len(vals) == 1:
                    val = vals[0]
                else:
                    val = tuple(vals)
            
            # parse source
            if source is not None:
                try:
                    source = sources[source]
                except KeyError:
                    raise ValueError("Table cell (%d, %d) has unknown source key: '%s'" % (i, j, source))
            
            cells[i][j] = (val, source)

    changes = [] # a list of parameters that are changed if we are rewriting a table
    for i,row in enumerate(row_names):
        for j,col in enumerate(col_names):
            kwds[row_key] = row
            kwds[col_key] = col
            oldval = setval(cells[i][j], name, **kwds)
            if oldval is not None and oldval != cells[i][j]:
                key = mk_key(name, **kwds)
                changes.append({'key': key, 'new': cells[i][j], 'old': oldval, 'name': name})
                #changes.append({'name': name, 'row': row, 'col': col, 'new': cells[i][j], 'old': oldval})
    return changes


def report_changes(changes):
    """
    For changes to data tables, give user a readout
    """
    if len(changes) > 0:
        print("\nWarning: Data Table '%s' (in memory) has been modified!" % changes[0]['name'])
        for ch in changes:
            # print('  >>> Changing %s, %s from default (%s) to %s' % (ch['row'], ch['col'], str(ch['new'][0]), str(ch['old'][0])))
            print('  >>> Changing %s, from default (%s) to %s' % (ch['key'], str(ch['old'][0]), str(ch['new'][0])))


def parse_sources(lines):
    sources = {}
    key = None
    val = []
    for l in lines:
        l = l.lstrip()
        m = re.match(r'\s*\[([^\]]+)\]\s+(.*)$', l)
        if m is not None:
            key = m.groups()[0]
            sources[key] = m.groups()[1].strip()
        else:
            if key is None:
                if l == '':
                    continue
                raise ValueError('Incorrect sources format--got text without '
                                 'citation index: "%s".' % l)
            sources[key] += '\n' + l
    return sources

#parse_sources('''\n\n[1] source 1\n    it's cool.\n[2] source 2 is not\n'''.split('\n'))




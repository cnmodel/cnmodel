from collections import OrderedDict
import re


# Unified collection point for all empirically-determined biophysical
# values. Each value is a tuple (val, source). 
DATA = OrderedDict()

def get(*args, **kwds):
    key = mk_key(*args, **kwds)
    return DATA[key][0]
    
def get_source(*args, **kwds):
    key = mk_key(*args, **kwds)
    return DATA[key][1]
    
def setval(key, val):
    DATA[key] = val
    
def mk_key(*args, **kwds):
    key = list(args) + list(kwds.items())
    key.sort(key=lambda a: a[0] if isinstance(a, tuple) else a)
    return tuple(key)


def add_table_data(name, row_key, col_key, data, **kwds):
    """
    Read data like
    
    '''
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
    '''
    
    
    """
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
    #
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
                raise Exception("Table line %d column %s does not obey column boundaries." % (i, j))
            
    # Break table into cells
    cells = []
    for line in table:
        if line.strip() != '':
            cells.append([line[cols[i]:cols[i+1]].strip() for i in range(len(cols)-1)])
    #print cells
    
    # Extract row/column names
    col_names = cells.pop(0)[1:]
    row_names = [cells[i].pop(0) for i in range(len(cells))]
    
    # Parse cell values
    for i in range(len(cells)):
        for j in range(len(cells[0])):
            cell = cells[i][j].strip()
            m = re.match(r'([^\[]*)(\[([^\]]+)\])?', cell)  # match like "0.7 [3]"
            if m is None:
                raise ValueError("Table cell (%d, %d) has bad format: '%s'" % (i, j, cell))
            
            # parse value
            val, _, source = m.groups()
            if val.strip() == '':
                val = None
            else:
                try:
                    val = int(val)
                except ValueError:
                    try:
                        val = float(val)
                    except ValueError:
                        raise ValueError("Table cell (%d, %d) value has bad format: '%s'" % (i, j, val))
            
            # parse source
            if source is not None:
                try:
                    source = sources[source]
                except KeyError:
                    raise ValueError("Table cell (%d, %d) has unknown source key: '%s'" % (i, j, source))
            
            cells[i][j] = (val, source)
    
    #print col_names
    #print row_names
    #print cells
    for i,row in enumerate(row_names):
        for j,col in enumerate(col_names):
            kwds[row_key] = row
            kwds[col_key] = col
            key = mk_key(name, **kwds)
            setval(key, cells[i][j])



def parse_sources(lines):
    sources = {}
    key = None
    val = []
    for l in lines:
        m = re.match(r'\s*\[([^\]])+\]\s+(.*)$', l)
        if m is not None:
            if key is not None:
                sources[key] = '\n'.join(val)
            elif len(val) > 0:
                raise ValueError('Incorrect sources format--got text without '
                                 'citation index "[N]": "%s".' % val)
            val = [m.groups()[1]]
            key = m.groups()[0]
        elif len(l.strip()) > 0:
            val.append(l)
    if key is not None:
        sources[key] = '\n'.join(val)
    elif len(val) > 0:
        raise ValueError('Incorrect sources format--got text without '
                         'citation index "[N]": "%s".' % val)
    return sources

#parse_sources('''\n\n[1] source 1\n    it's cool.\n[2] source 2 is not\n'''.split('\n'))




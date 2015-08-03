# -*- encoding: utf-8 -*-
import pytest
from cnmodel import data


table = u"""

Description of data
It has multiple lines

And empty lines.

-------------------------------------------------------------------------------
          col1         col2           col3
param1    15±6.5 [1]   2.2±1.5 [2]    0.87±0.23 [3]
param2    3            5              0.87±0.23 [3]
param3      3.4                       7 [2]
param4    1 [12]           
-------------------------------------------------------------------------------

[1] Source 1
[2] Multiline source
        #2
    end of #2
[3] Multiline source
        #3
    end of #3

[12] another
     source


"""

data.add_table_data('test_data', row_key='param', col_key='col', data=table,
                    extra='test_kwd')


def test_db():
    with pytest.raises(TypeError):
        # raise exception if unicode is given in ono-unicode string
        data.add_table_data('test_data', row_key='param', col_key='col', data='±')
    
    d = data.get('test_data', param='param1', col='col2', extra='test_kwd')
    assert d == (2.2, 1.5)

    d = data.get('test_data', param='param1', col=['col1', 'col2', 'col3'], extra='test_kwd')
    assert d == {'col1': (15, 6.5), 'col2': (2.2, 1.5), 'col3': (0.87, 0.23)}
    
    d = data.get('test_data', param='param2', col=['col1', 'col2', 'col3'], extra='test_kwd')
    assert d == {'col1': 3, 'col2': 5, 'col3': (0.87, 0.23)}
    
    d = data.get('test_data', param='param3', col=['col1', 'col2', 'col3'], extra='test_kwd')
    assert d == {'col1': 3.4, 'col2': None, 'col3': 7}
    
    s = data.get_source('test_data', param='param1', col='col2', extra='test_kwd')
    assert 'end of #2' in s
    
    s = data.get_source('test_data', param='param2', col='col2', extra='test_kwd')
    assert s is None

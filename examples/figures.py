from __future__ import print_function
"""

"""
import sys
import subprocess

if len(sys.argv) < 2:  # if no argument, print helpful message
    print("Plot selected figures from paper, Manis and Campagnoal, Hearing Research. 2018")
    print("Usage: figures.py [2a | 2b | 2c | 3 | 4 | 5 | 6a | 6d | 7]")
    exit(1)
    
arg = sys.argv[1]  # get argument, check that it is valid
if arg not in ['2a', '2b', '2c', '3', '4', '5', '6a', '6d', '7']:
    print("Usage: figures.py [2a | 2b | 2c | 3 | 4 | 5 | 6a | 6d | 7]")
    exit(1)
    
if arg == '2a':
    proc = subprocess.Popen(['python', 'examples/test_mechanisms.py', 'klt'],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
    proc.wait()
    print (proc.stdout.read())
     #    ;;
     # 2b)
     #    python examples/test_mechanisms.py kht
     #    ;;
     # 2c)
     #    python examples/test_mechanisms.py ka
     #    ;;
     # 3)
     #     python examples/toy_model.py
     #     ;;
     # 4)
     #     python examples/test_synapses.py sgc bushy
     #     ;;
     # 5)
     #     python examples/test_decorator.py
     #     ;;
     # 6a)
     #     python examples/test_bushy_variation.py a
     #     ;;
     # 6d)
     #     python examples/test_bushy_variation.py d
     #     ;;
     #
     # 7)
     #     while true; do
     #         echo "This figure may take hours to generate!"
     #         read -p "Are you sure you want to run the script?" yn
     #         case $yn in
     #             [Yy]* ) python examples/test_physiology.py; break;;
     #             [Nn]* ) exit;;
     #             * ) echo "Please answer yes or no.";;
     #         esac
     #     done
     #     ;;
     #



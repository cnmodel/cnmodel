case $1 in
     2a)
        python examples/test_mechanisms.py klt
        ;;
     2b)
        python examples/test_mechanisms.py kht
        ;;
     2c)
        python examples/test_mechanisms.py ka
        ;;
     3)
         python examples/toy_model.py
         ;;
     4) 
         python examples/test_synapses.py sgc bushy
         ;;
     5)
         python examples/test_decorator.py
         ;;
     6a)
         python examples/test_bushy_variation.py a
         ;;
     6b)
         python examples/test_bushy_variation.py b
         ;;

     7)
         while true; do
             echo "This figure may take hours to generate!"
             read -p "Are you sure you want to run the script?" yn
             case $yn in
                 [Yy]* ) python examples/test_physiology.py; break;;
                 [Nn]* ) exit;;
                 * ) echo "Please answer yes or no.";;
             esac
         done
         ;;
         
     *)
         echo $"Plot selected figures from paper"
         echo $"Usage: $0 {2a | 2b | 2c | 3 | 4 | 5 | 6a | 6b | 7}"
         exit 1         
esac


FLICK_FOLDER=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
if [ -n "$($SHELL -c 'echo $ZSH_VERSION')" ]; then
    echo 'export FLICK_PATH='$FLICK_FOLDER >> ~/.zprofile
    echo 'export FLICK_COMPILER="clang++ -std=c++20"' >> ~/.zprofile
    echo 'export PATH=$FLICK_PATH/main:$PATH' >> ~/.zprofile
    echo 'Three lines have been added to your ~/.zprofile'
    echo 'Now, run the command'
    echo ''
    echo '  source ~/.zprofile'
    echo ''
    echo 'and then rerun make'
    echo ''
    exit 1
    
elif [ -n "$($SHELL -c 'echo $BASH_VERSION')" ]; then
    echo 'export FLICK_PATH='$FLICK_FOLDER >> ~/.bashrc
    echo 'export FLICK_COMPILER="g++ -std=c++20"' >> ~/.bashrc
    echo 'export PATH=$FLICK_PATH/main:$PATH' >> ~/.bashrc
    echo 'Three lines have been added to your ~/.bashrc'
    echo 'Now, run the command'
    echo ''
    echo '  source ~/.bashrc'
    echo ''
    echo 'and then rerun make'
    echo ''
    exit 1

else
    echo 'Unknown shell. Set these environmental variables in your'
    echo 'shell profile script, source it, and run make again'
    echo ''
    echo 'FLICK_PATH=$HOME/Flick-RT'
    echo 'FLICK_COMPILER="g++ -std=c++20"'
    echo 'PATH=$FLICK_PATH/main:$PATH'
    echo ''
    exit 1
fi



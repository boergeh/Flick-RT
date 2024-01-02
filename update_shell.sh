FLICK_FOLDER=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
if [ -n "$($SHELL -c 'echo $ZSH_VERSION')" ]; then
    echo 'export FLICK_PATH='$FLICK_FOLDER >> ~/.zprofile
    echo 'export FLICK_COMPILER="clang++ -std=c++20"' >> ~/.zprofile
    echo 'export PATH=$FLICK_PATH/main:$PATH' >> ~/.zprofile
    echo 'Three lines have been added to your ~/.zprofile. Now, rerun make.'
    exec zsh
    
elif [ -n "$($SHELL -c 'echo $BASH_VERSION')" ]; then
    echo 'export FLICK_PATH='$FLICK_FOLDER >> ~/.profile
    echo 'export FLICK_COMPILER="g++ -std=c++20"' >> ~/.profile
    echo 'export PATH=$FLICK_PATH/main:$PATH' >> ~/.profile
    echo 'Three lines have been added to your ~/.profile. Now, rerun make.'
    exec bash

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



export FLICK_FOLDER=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
if [ -n "$($SHELL -c 'echo $ZSH_VERSION')" ]; then
    echo 'Need to add flick environmental variables to your ~/.zprofile'
    read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY]
	|| $confirm == [yY][eE][sS] ]] || exit 1
    echo 'export FLICK_PATH='$FLICK_FOLDER >> ~/.zprofile
    echo 'export FLICK_COMPILER="clang++ -std=c++20"' >> ~/.zprofile
    echo 'export PATH=$FLICK_PATH/main:$PATH' >> ~/.zprofile
    #zsh
    #source  ~/.zprofile
    #$($SHELL -c '~/.zprofile')
    
elif [ -n "$($SHELL -c 'echo $BASH_VERSION')" ]; then
    echo 'Need to add flick environmental variables to your ~/.profile'
    read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY]
	|| $confirm == [yY][eE][sS] ]] || exit 1
    echo 'export FLICK_PATH='$FLICK_FOLDER >> ~/.profile
    echo 'export FLICK_COMPILER="g++ -std=c++20"' >> ~/.profile
    echo 'export PATH=$FLICK_PATH/main:$PATH' >> ~/.profile
    source  ~/.profile
else
    echo 'Unknown shell!'
    echo 'Set these environmental variables in your shell profile'
    echo 'and run make again'
    echo ''
    echo 'FLICK_PATH=$HOME/Flick-RT'
    echo 'FLICK_COMPILER="clang++ -std=c++20"'
    echo 'PATH=$FLICK_PATH/main:$PATH'
    echo ''
    exit 1
fi



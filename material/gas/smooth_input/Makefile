SRC      = $(wildcard *.cpp)
OBJ      = $(SRC:.cpp=.o)
DEP      = $(patsubst %.cpp,%.d,$(SRC))
NAME     = $(SRC:.cpp=)

all:	obj link test

clean:
	@rm -f $(OBJ) $(DEP) ./$(NAME) ./*~ ./obj ./tmp.txt
link:
	$(FLICK_COMPILER) -o $(NAME) $(OBJ)

obj:	$(OBJ)

test:
	@./$(NAME)

-include $(DEP) 

%.o: %.cpp Makefile
	$(FLICK_COMPILER) -MMD -MP -c $< -o $@

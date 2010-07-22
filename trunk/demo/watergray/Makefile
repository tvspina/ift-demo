# Diretórios de acesso 

INC=../../include
LIB=../../lib
BIN=.

INCLUDES= -I$(INC) -I$(OPF_DIR)/include

#FLAGS= -g -Wall 
FLAGS= -O3 -Wall

$(LIB)/libift.a:
	cd ../../; make libift; cd -


 $@.c: $(LIB)/libift.a $@.c 
	 gcc $(FLAGS) $@.c -o $(BIN)/$@ $(INCLUDES) -L$(LIB) -lift -lm ; 












#Directories
BIN_DIR  = ./bin
EXEC_DIR = ./src/backend/lib
OBJ_DIR  = ./obj
SRC_DIR  = ./src/backend/cpp

#Compiler
CXX = g++
CFLAGS = -g -std=c++11
INCLUDES = -I $(SRC_DIR) -I /usr/include/python3.5m -I ./eigen-eigen-323c052e1731/
LFLAGS = -lboost_python-py35 -lpython3.5m -ldl -lppl -lglpk -lgmp

CPP_SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
CPP_OBJECTS := $(CPP_SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

#Target
EXEC     = libc2e2.so
TARGET  :=  $(EXEC_DIR)/$(EXEC)

.PHONY: all
all: $(TARGET)

.PHONY: dirs
dirs:
	@ mkdir -p $(OBJ_DIR)
	@ mkdir -p $(BIN_DIR)

# FIXME
$(TARGET): $(CPP_OBJECTS)
	$(CXX) -shared -Wl,--no-undefined $(CFLAGS) -o $@ $(CPP_OBJECTS) $(LFLAGS) $(INCLUDES)

# FIXME
$(CPP_OBJECTS): $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CFLAGS) -fPIC -c -o $@ $< $(INCLUDES)
	#$(CXX) $(CFLAGS) $(INCLUDES) -fPIC -c -o $@ $< 

.PHONY: clean
clean:
	rm -rf $(OBJ_DIR)/*.o
	rm -rf $(EXEC_DIR)/*.so

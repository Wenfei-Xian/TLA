CXX = g++

CXXFLAGS = -Wall

TARGET = TLA

SRC = TLA.cpp

INCLUDES = -I$(HTSLIB_PATH)/include \
           -I$(SPOA_PATH)/include

LDFLAGS = -L$(HTSLIB_PATH)/lib \
          -L$(SPOA_PATH)/build/lib

LIBS = -lhts -lspoa

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(INCLUDES) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(TARGET)

.PHONY: all clean

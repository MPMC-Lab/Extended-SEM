FC = ifort

TARGET = a.out

BUILD = global.o main.o

#global = global.o
#main = main.o
#OBJS = $(global) $(main)

$(TARGET): $(BUILD)
	$(FC) $(BUILD) -o $(TARGET)

%.o: %.f90
	$(FC) -c $<

all : $(TARGET)

clean :
	rm *.mod *.o *.out

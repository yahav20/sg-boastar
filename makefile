CFLAGS      = -Wall -O3
LFLAGS      = -Wall -O3
LDLIBS      = -lm
CC      = gcc

OBJ     = obj/heap.o obj/boastar.o obj/graph.o obj/main_boa.o
OBJ_BOD = obj/heap.o obj/bod.o obj/graph.o obj/main_bod.o
OBJ_NAMOADR = obj/heap.o obj/namoadr.o obj/graph.o obj/main_namoadr.o
OBJ_SG_BOA = obj/heap.o obj/boastar.o obj/sg_boastar.o obj/graph.o obj/main_sg_boa.o
OBJ_BENCHMARK_SG_BOA = obj/heap.o obj/boastar.o obj/sg_boastar.o obj/graph.o obj/benchmark_sg_boa.o

all: boa bod namoadr sg_boa benchmark_sg_boa

boa:  $(OBJ)
	$(CC) $(LFLAGS) -o boa $(OBJ)

bod:  $(OBJ_BOD)
	$(CC) $(LFLAGS) -o bod $(OBJ_BOD)

namoadr:  $(OBJ_NAMOADR)
	$(CC) $(LFLAGS) -o namoadr $(OBJ_NAMOADR)

sg_boa: $(OBJ_SG_BOA)
	$(CC) $(LFLAGS) -o sg_boa $(OBJ_SG_BOA) $(LDLIBS)

benchmark_sg_boa: $(OBJ_BENCHMARK_SG_BOA)
	$(CC) $(LFLAGS) -o benchmark_sg_boa $(OBJ_BENCHMARK_SG_BOA) $(LDLIBS)

obj/%.o: src/%.c
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -f obj/*.o boa bod namoadr sg_boa benchmark_sg_boa

rebuild: clean all

test-all: test-NY test-BAY test-NE

test-NY:
	./benchmark_sg_boa Maps/NY-road-d.txt 20 60	
	
test-BAY:
	./benchmark_sg_boa Maps/BAY-road-d.txt 20 60

test-NE:
	./benchmark_sg_boa Maps/NE-road-d.txt 20 60

test-all: test-NY test-BAY test-NE
CC = cc
CFLAGS = -lm -lpthread -g

all:
	$(CC) -o mtransport $(CFLAGS) src/main.c src/transport.c src/reactions.c src/tmath.c src/input.c

clean:
	rm mtransport

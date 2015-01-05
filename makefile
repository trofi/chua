TARGET := chua_gen
CFLAGS := -Wall -Wextra -O2
LIBS   := -lglut -lGL -lm

SOURCES := $(wildcard *.c)
OBJECTS := $(SOURCES:.c=.o)

all: $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@

.PNONY: clean run

clean:
	@rm -rf $(TARGET) $(OBJECTS) *~ semantic.cache *TAGS *BROWSE

run: $(TARGET)
	./$(TARGET)

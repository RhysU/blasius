CFLAGS   := -std=c99 -Wall -Wextra -O3
CPPFLAGS := `gsl-config --cflags`
LDFLAGS  := `gsl-config --libs`

all     : blasius
blasius : blasius.c

clean:
	rm -f blasius

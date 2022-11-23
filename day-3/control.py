#!/usr/bin/env python

from time import sleep
from math import log2

command = 'whatever'


while True:
    command = input('what is your command? ')
    if command == 'quit':
        break
    elif command == 'sleep':
        sleep(3)
    elif command == 'eat':
        print('yum')

    elif command.startswith('add '):
        try:
            number1, number2 = command.split()[1:]
        except ValueError:
            print("Hey, you fool, I can only add two numbers! Try again.")
            continue

        print(int(number1) + int(number2))

    elif command.startswith('log '):
        _, number1 = command.split()
        print(log2(float(number1)))
    else:
        print(f'You said {command!r}')

print('Goodbye!')

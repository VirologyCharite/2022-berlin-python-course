#!/usr/bin/env python


pegs = [
    [8, 7, 6, 5, 4, 3, 2, 1],
    [],
    [],
]

moves = 0


def drawPegs(move):
    print(f'After move {move}:')
    for pegNumber, peg in enumerate(pegs):
        print('Peg', pegNumber)
        if peg:
            for disc in reversed(peg):
                print('  ', '-' * disc, f'({disc})')
        else:
            print('  empty')

    print()


def move(origin, destination, nDiscs):
    """
    Move all the pegs on the origin peg to the destination peg.
    """
    global moves
    moves += 1

    drawPegs(moves)


drawPegs(0)
move(0, 1, 8)

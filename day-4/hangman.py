import sys
import os
import getpass

def generateHangy(guess_count, max_guesses):
    step = int(guess_count/max_guesses*4)

    if step < 2:
        hangy = """
|
|
|
        
        """

    elif step < 3:
        hangy = """
|
|
|
|
|
|
        
        """
    elif step < 4:
        hangy = """
 ___
|   
|   
|  
|  
|  
|
        """
    elif step == 5:
        hangy = """
 ___
|   |
|   O
|  /|\\
|   
|  
|
        """  
    elif step == 6:
        hangy = """
 ___
|   |
|   O
|  /|\\
|   |
|  / \\
|
        """  
    return(hangy)   


secret_word = getpass.getpass("secret word: ").lower()

max_guesses = int(input("How many guesses should be allowed? "))

guessed = set()


os.system('clear')



print(f"The word has {len(secret_word)} characters")


guess_count = 0
while True:
    guess = input("Guess a letter or word: ").lower()
    guess_count += 1
    if len(guess) > 1:
        if (guess == secret_word):
            print("Yey you got it")
            revealed_count = len(secret_word)
            break 
    guessed.add(guess)
    revealed_word = ""
    revealed_count = 0
    for letter in secret_word:
        if letter in guessed:
            revealed_word += letter
            revealed_count += 1
        else:
            revealed_word += "_"
    print(revealed_word)
    print(f"revealed {revealed_count} of {len(secret_word)}")
    print(f"guesses left {max_guesses-guess_count}")
    print(generateHangy(guess_count, max_guesses))
    if revealed_count == len(secret_word):
        print("you are a winner")
        print()
        break
    elif guess_count == max_guesses:
        print("you are a looser")
        break



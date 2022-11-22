def celciusToFahrenheit(c):
    return (9 / 5) * c + 32


def say(mesg):
    print("The message:", mesg)


for temp in range(-45, -10):
    say(temp)
    print(temp, "C = ", celciusToFahrenheit(temp), "F")

# for number1 in 1, 2, 3, 4:

for number1 in list(range(1, 4)):
    print("First number is now", number1)

    for number2 in 1, 2, 3:
        print("  Second number is now", number2)

        for number3 in 1, 2:
            print("    ", number1, "x", number2, "x", number3, "=",
                  number1 * number2 * number3)

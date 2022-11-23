import sys

wordCounts = {}

for line in sys.stdin:
    for word in line.split():
        if word in wordCounts:
            # This word is already in the dictionary. Add one to its count.
            wordCounts[word] = wordCounts[word] + 1
        else:
            # This is a new word, not already in the dictionary.
            wordCounts[word] = 1

maxCount = 0
mostCommonWords = []

for word in wordCounts:
    count = wordCounts[word]
    if count > maxCount:
        maxCount = count
        mostCommonWords = [word]
    elif count == maxCount:
        mostCommonWords.append(word)


if len(mostCommonWords) > 1:
    # There are multiple words that are most common.
    words = ", ".join(sorted(mostCommonWords))
    print(f"There are {len(mostCommonWords)} words that are most common: "
          f"{words} with count {maxCount}.")
else:
    # We know there is only one word.
    word = mostCommonWords[0]
    print(f"The most common word is {word!r} with count {maxCount}.")

# Equivalently:
# print("The most common word is", mostCommonWord, "with count", maxCount)

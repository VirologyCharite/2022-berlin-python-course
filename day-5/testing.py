#!/usr/bin/env python


def orig_cleanWord(word):

    for punc in ",$-()!":
        word = word.replace(punc, "")

    return word


def cleanWord(word):

    clean = ""

    for letter in word:
        if letter.isalpha():
            clean += letter

    return clean


def testEmptyString():
    "If we pass the empty string we should get the empty string back."
    assert cleanWord("") == ""


def testCheckWordWithNoPunctuation():
    assert cleanWord("hello") == "hello"


def testRemoveCommaAtEnd():
    assert cleanWord("hello,") == "hello"


def testRemoveCommasAtEnd():
    assert cleanWord("hello,,,") == "hello"


def testRemoveDollarSignsAtEnd():
    assert cleanWord("hello$$") == "hello"

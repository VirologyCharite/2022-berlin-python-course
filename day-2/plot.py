import matplotlib.pyplot as plt

coffee = [
    6, 3, 4, 1, 25, 11
]

iq = [
    100, 45, 55, 20, 130, 50
]


figure, ax = plt.subplots()

# ax.scatter(coffee, iq)
# ax.plot(coffee, iq, linewidth=1, color="red", marker="o", markersize=10)
ax.plot(coffee, iq, "bo")

ax.set_xlabel("Coffee (cups per day)")
ax.set_ylabel("Intelligence quotient")
ax.set_title("The title of the scatter plot")


figure.suptitle("Here is the title")
# figure.legend()

plt.savefig("figure.pdf")

# plt.show()

class Bag:
    "This is a class for bags."

    def __init__(self, waterproof=False, brand='Dakine',
                 height=10, depth=10, width=10):
        self.waterproof = waterproof
        self.brand = brand
        self.height = height
        self.depth = depth
        self.width = width

    def __str__(self):
        if self.waterproof:
            return 'Congratulations on buying a waterproof %s bag' % self.brand
        else:
            return 'Congratulations on buying a %s bag' % self.brand

    def volume(self):
        return self.height * self.depth * self.width

###


b1 = Bag(waterproof=True)
print(b1)

b2 = Bag(waterproof=False)
print(b2)

b3 = Bag(height=3)
print(b3)

print('Bag3 has volume', b3.volume())

"""
48 channel bolometer system, in 4 arrays (A,B,C,D)

"""

# Definitions of the bolometer arrays, their apertures and channels
arrays = {
    "A": {
        "Aper": [0.6985, 0.70104],
        # Locations are [X (cm), Y (cm), Z (cm), angle (degrees)]
        "Aper Loc": [234.9529342, 72.99039831, 2.584120859, 236.2962527],
        "channels": {"#1": [235.0679877, 83.90805959, 2.12571092, 269.3249741]},
    }
}

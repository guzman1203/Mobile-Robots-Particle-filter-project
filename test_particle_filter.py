import pytest
from particle_filter_helper import (
    get_environment_coords, 
    get_cardinal_between_neighbors, 
    get_all_shortest_paths,
    get_shortest_path
)

@pytest.mark.parametrize(
    "cell_number, expected",
    [
        ( 1,(2,-2)), # Cell 1 at    2, -2
        (13,(0,0)), # Cell 13 at   0,  0
        (25,(-2,2)) # Cell 25 at  -2,  2
    ]
)
def test_get_coordinates_from_cellnumber(cell_number, expected):
    assert get_environment_coords(cell_number) == expected

@pytest.mark.parametrize(
    "start, neighbor, expected",
    [
        ( 1,  2, "E"), # East 
        (13,  8, "N"), # North
        (25, 24, "W"), # West
        (14, 19, "S")  # South
    ]
)
def test_get_cardinal_between_neighbors(start, neighbor, expected):
    assert get_cardinal_between_neighbors(start, neighbor) == expected

def test_get_all_shortest_paths():
    start = 13
    expected = {
         1: (4, "WWNN"),  2: ( 5, "WWNNE"),         3: ( 6, "WWNNEE"),       4: ( 7, "WWNNEEE"),     5: (8, "WWNNEEEE"),
         6: (3, "WWN"),   7: (12, "WWNNEEEESWWW"),  8: (11, "WWNNEEEESWW"),  9: (10, "WWNNEEEESW"), 10: (9, "WWNNEEEES"),
        11: (2, "WW"),   12: ( 1, "W"),            13: ( 0, ""),            14: ( 1, "E"),          15: (2, "EE"),
        16: (3, "WWS"),  17: ( 4, "WWSE"),         18: ( 5, "WWSEE"),       19: ( 4, "EESW"),       20: (3, "EES"),
        21: (4, "WWSS"), 22: ( 5, "WWSSE"),        23: ( 6, "WWSSEE"),      24: ( 5, "EESSW"),      25: (4, "EESS"), 
    }
    assert get_all_shortest_paths(start) == expected

def test_get_shortest_path():
    start, goal = 7, 21
    expected = "EEENWWWWSSSS"
    assert get_shortest_path(start, goal) == expected

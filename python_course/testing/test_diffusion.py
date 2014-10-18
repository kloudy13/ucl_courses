
from nose.tools import assert_true, assert_equal
from diffusion_model import *

def test_my_function_exact_input():
    particles = [1,3,7,2,9];
    D = 2.0;
    result = energy(particles, D);
    assert_equal(result, 122.0);



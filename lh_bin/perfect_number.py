from iga.apps.base import emain
import logging


def get_divisors(n):
    """
    n: the input number to get divisors
    """
    for i in (range(1, n//2 + 2)):
        if n % i == 0:
            yield i


def is_perfect_num(n=None):
    n = int(n)
    sum = 0
    divisor_list = []
    for divisor in get_divisors(n):
        #logging.debug(divisor)
        sum += divisor
        divisor_list.append(divisor)
    if sum == n:
        return divisor_list
    else:
        return 0


def find_perfect_num(min=1, max=1000):
    for i in range(min, max + 1):
        is_perfect = is_perfect_num(i)
        if is_perfect:
            logging.info("Found pefect number {} with divisors {}".format(i, is_perfect))
            print(i)


if __name__ == "__main__":
    emain()

GREY = '\033[39m'
RED = '\033[31m'
GREEN = '\033[34m'


# TESTS
class Test:

    def describe(name):
        print("SET:", name)

    def it(name):
        print("Testing", name)

    def assert_equals(self, a, b, error_message=""):
        if a == b:
            print(GREEN, "Okay!", a, "equals", b)
        else:
            print(RED, "Bad!", a, "!=", b, ".", error_message)
        print(GREY)

    def expect_error(self, *kwargs):
        pass

test = Test()

from PQAnalysis.utils.decorators import count_decorator, instance_function_count_decorator


def test_count_decorator():
    @count_decorator
    def test_func(reset_counter=False):
        pass

    assert test_func.counter == 0
    test_func()
    assert test_func.counter == 1
    test_func()
    assert test_func.counter == 2
    test_func(reset_counter=True)
    assert test_func.counter == 1
    test_func()
    assert test_func.counter == 2


def test_instance_function_count_decorator():
    class Class:
        def __init__(self):
            self.instance_counter = 0

        @instance_function_count_decorator
        def test_func(self, reset_counter=False):
            self.instance_counter = self.counter[Class.test_func.__name__]

    instance = Class()
    assert instance.instance_counter == 0
    instance.test_func()
    assert instance.instance_counter == 1
    instance.test_func()
    assert instance.instance_counter == 2
    instance.test_func(reset_counter=True)
    assert instance.instance_counter == 1
    instance.test_func()
    assert instance.instance_counter == 2

    instance2 = Class()
    assert instance2.instance_counter == 0
    instance2.test_func()
    assert instance2.instance_counter == 1
    instance2.test_func()
    assert instance2.instance_counter == 2

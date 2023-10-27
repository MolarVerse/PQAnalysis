from PQAnalysis.utils.decorators import count_decorator


def test_count_decorator():
    @count_decorator
    def test_func():
        pass

    assert test_func.counter == 0
    test_func()
    assert test_func.counter == 1
    test_func()
    assert test_func.counter == 2
    test_func()
    assert test_func.counter == 3

from beartype.typing import Any
from PQAnalysis.formats import BaseEnumFormat



def test_BaseEnumFormat():

    class TestEnum(BaseEnumFormat):
        A = "a"
        B = "b"

        @classmethod
        def _missing_(cls, value: object) -> Any:
            return super()._missing_(value, ValueError)

    assert TestEnum.A == "a"
    assert TestEnum.B == "b"
    assert TestEnum("a") == TestEnum.A
    assert TestEnum("b") == TestEnum.B
    assert TestEnum("b") != 2

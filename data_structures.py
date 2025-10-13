from typing import Self


class IntPair:

    def __init__(
        self,
        elt_a: int,
        elt_b: int
    ) -> Self:
        self.start = elt_a if elt_b >= elt_a else elt_b
        self.end = elt_b if elt_b >= elt_a else elt_a

    def __str__(self) -> str:
        return f'[{self.start}-{self.end}]'

    def __repr__(self) -> str:
        return f'[{self.start}-{self.end}]'


class Interval:
    """
    Class that implements basic operations on intervals in python.

    It offers the following:

    **Overriding common operators.**
    + The plus sign will stand for union
    + The mul sign for intersection
    + The minus sign for symmetric difference
    """

    def __init__(
        self,
        intervals: list[IntPair] | IntPair,
    ) -> Self:
        self.__intervals = [intervals] if isinstance(
            intervals, IntPair) else sorted(intervals, key=lambda x: x.start)
        # Computing union of intervals
        new_intervals: list[IntPair] = list()
        i: int = 0
        found_ending: bool = True
        while i < len(self.__intervals):
            if found_ending:
                start, end = self.__intervals[i].start, self.__intervals[i].end

            if self.__intervals[i].start <= end and start <= self.__intervals[i].end:
                start = min(self.__intervals[i].start, start)
                end = max(self.__intervals[i].end, end)
                i += 1
                found_ending: bool = False
            else:
                new_intervals.append(IntPair(start, end))
                i += 1
                found_ending: bool = True
        if not found_ending:
            new_intervals.append(IntPair(start, end))
        self.intervals = sorted(new_intervals, key=lambda x: x.start)
        del self.__intervals

    def total(self) -> int:
        """Returns the total count values of Interval

        Returns:
            int: number of values in interval
        """
        return sum([i.end-i.start for i in self.intervals])

    def symmetric_difference(self, elt: Self) -> Self:
        raise NotImplementedError()

    def intersect(self, elt: Self) -> Self:
        new_intervals: list[IntPair] = list()
        #start, end = float('-inf'), float('inf')
        i, j = 0, 0
        while i < len(self.intervals)-1 and j < len(elt.intervals)-1:
            if self.intervals[i].start <= elt.intervals[j].end and elt.intervals[j].start <= self.intervals[i].end:
                new_intervals.append(IntPair(max(self.intervals[i].start, elt.intervals[j].start), min(self.intervals[i].end, elt.intervals[j].end)))
                i+=1
                j+=1
            else:
                if self.intervals[i].end < self.intervals[j].end:
                    i+=1
                else:
                    j+=1
        return Interval(new_intervals)

    def union(self, elt: Self) -> Self:
        new_intervals: list[IntPair] = list()
        start, end = float('inf'), float('-inf')
        i, j = 0, 0
        while i < len(self.intervals) and j < len(elt.intervals):
            if self.intervals[i].start <= elt.intervals[j].end and elt.intervals[j].start <= self.intervals[i].end:
                start = min(self.intervals[i].start, elt.intervals[j].start)
                end = max(self.intervals[i].end, elt.intervals[j].end)
            else:
                new_intervals.append(IntPair(start, end))
                start, end = float('inf'), float('-inf')
        return Interval(new_intervals)
    
    def to_set(self) -> set:
        return {x for intpair in self.intervals for x in range(intpair.start,intpair.end)}

    def __add__(self, elt: Self) -> Self:
        return self.union(elt)

    def __radd__(self, elt: Self) -> Self:
        return self.union(elt)

    def __sub__(self, elt: Self) -> Self:
        return self.symmetric_difference(elt)

    def __rsub__(self, elt: Self) -> Self:
        return self.symmetric_difference(elt)

    def __mul__(self, elt: Self) -> Self:
        return self.intersect(elt)

    def __rmul__(self, elt: Self) -> Self:
        return self.intersect(elt)

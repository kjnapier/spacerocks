from typing import List

class Detection:

    def __init__(self):
        pass


class Link:

    def __init__(self, detections: List[Detection]):
        self.detections = detections


class OrbitFitter:

    def __init__(self, links: Link):
        self.links = links

    def fit(self):
        pass
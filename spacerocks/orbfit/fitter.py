from ..linking.detection import Detection


class Link:

    def __init__(self, detections: list[Detection]):
        self.detections = detections


class OrbitFitter:

    def __init__(self, links: Link):
        self.links = links

    def fit(self):
        pass
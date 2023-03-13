from .detection import Detection


class Link:

    def __init__(self, detections: list[Detection]):
        self.detections = detections
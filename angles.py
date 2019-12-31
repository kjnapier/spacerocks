# take an angle as input, and return in the desired form

def DegreesToHMS(angle):
    a = angle / 15
    h = int(a)
    m = int(60 * (a - h))
    s = (60 * (a - h) - m) * 60
    return '{}h {}m {}s'.format(h, m, s)

def DegreesToDMS(angle):
    d = int(angle)
    m = int(60 * (angle - d))
    s = (60 * (angle - d) - m) * 60
    return '{}:{}:{}'.format(d, m, s)
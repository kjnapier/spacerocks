# Take an input date, and return the date in the desired form.

def JD(day, month, year, UT):
    
    if month > 2:
        y = year
        m = month
    else:
        y = year - 1
        m = month + 12
        
    if year > 1582:
        B = int(y / 400) - int(y / 100)
    elif year < 1582:
        B = -2
    else:
        if month < 10:
            B = -2
        elif month > 10:
            B = int(y / 400) - int(y /100)
        else:
            if day <= 4:
                B = -2
            elif day >= 15:
                B = int(y / 400) - int(y / 100)
                
    return int(365.25 * y) + int(30.6001 * (m + 1)) + 1720996.5 + B + day + UT / 24;

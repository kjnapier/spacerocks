# The `Units` Object

In my experience, units are one of the most common sources of bugs in 
solar system computations (or any computations, really). `spacerocks`
takes proactive measures to avoid such bugs by using a `Units` object.

Using a `Units` object allows you to specify a set of units once, localizing
any typos or mistakes to a single instance. You can then pass your `Units`
object to whatever class instantiation or method you are using.

```Python
from spacerocks.units import Units

units = Units()
print(units.current())
```

```zsh
>>>  Quantity             Unit           
     ---------------------------------------
     distance             AU             
     angle                deg            
     timescale            utc            
     timeformat           None           
     speed                AU / d         
     size                 km             
     density              g / cm3        
     mass                 kg             
     ra                   deg            
     dec                  deg 
```

You can then set the attributes of the `Units` object like this.

```Python
from astropy import units as u

units.timescale = 'tdb'
units.timeformat = 'jd'
units.angle = u.rad
```

While most fields can be set using raw strings (i.e. `'deg'` or `'rad'`), 
I recommend explicitly passing an `astropy` unit where applicable. 
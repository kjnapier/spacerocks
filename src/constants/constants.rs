use lazy_static::lazy_static;
use std::collections::HashMap;
use nalgebra::Matrix3;

pub const KM_TO_AU: f64 = 1.0 / 149_597_870.700;
pub const M_TO_AU: f64 = KM_TO_AU / 1000.0;

pub const SECONDS_PER_DAY: f64 = 86_400.0;

pub const EQUAT_RAD: f64 = 6378137.0;
pub const FLATTEN: f64 = 1.0 / 298.257223563;
pub const O_M_FLATTEN: f64 = 1.0 - FLATTEN;
pub const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

pub const MU_BARY: f64 = 0.00029630927493457475;
pub const SPEED_OF_LIGHT: f64 = 173.14463268466926; // speed of light in au/day

pub const GRAVITATIONAL_CONSTANT: f64 = 0.00029591220828559104;
// pub const GRAVITATIONAL_CONSTANT: f64 = 0.00029591220819207774;

pub const ROTATION_J2000: Matrix3<f64> = Matrix3::new(1.0, 0.0, 0.0,
                                                      0.0, 1.0, 0.0,
                                                      0.0, 0.0, 1.0);

pub const ROTATION_FK4: Matrix3<f64> = Matrix3::new(0.99992567949568767, 0.011181483239171792, 0.0048590037723143858,
                                                   -0.01118148322046629, 0.99993748489331347, -2.7170293744002029e-05,
                                                   -0.0048590038153592712, -2.7162594714247048e-05, 0.9999881946023742);

pub const ROTATION_GALACTIC: Matrix3<f64> = Matrix3::new(-0.054875539395742523, -0.87343710472759606, -0.48383499177002515,
                                                         0.49410945362774389, -0.44482959429757496, 0.74698224869989183,
                                                         -0.86766613568337381, -0.19807638961301985, 0.45598379452141991);

pub const ROTATION_ECLIPJ2000: Matrix3<f64> = Matrix3::new(1.0, 0.0, 0.0,
                                                           0.0, 0.91748206206918181, 0.39777715593191371,
                                                           0.0, -0.39777715593191371, 0.91748206206918181);

pub const ROTATION_INVARIABLE: Matrix3<f64> = Matrix3::new(-0.3023595432982142,  0.8743824933349968,  0.3795182677654754,
                                                           -0.9527924404431956, -0.2888105590660314, -0.0936832539501452,
                                                            0.0263975536104876, -0.389928162416098 ,  0.920429658444365);

// make a hash map of the rotation matrices
lazy_static! {
    pub static ref ROTATION_MATRICES: HashMap<String, Matrix3<f64>> = {
        let mut m = HashMap::new();
        m.insert("J2000".to_string(), ROTATION_J2000);
        m.insert("FK4".to_string(), ROTATION_FK4);
        m.insert("GALACTIC".to_string(), ROTATION_GALACTIC);
        m.insert("ECLIPJ2000".to_string(), ROTATION_ECLIPJ2000);
        m.insert("INVARIABLE".to_string(), ROTATION_INVARIABLE);
        m
    };
}

const KM_PER_AU: f64 = 149597870.700;
// const SECONDS_PER_DAY: f64 = 86400.0
const KM3_PER_SECOND2_TO_AU3_PER_DAY2: f64 = (1.0 / KM_PER_AU) * (1.0 / KM_PER_AU) * (1.0 / KM_PER_AU) * (SECONDS_PER_DAY * SECONDS_PER_DAY);

// lazy_static! {
//     pub static ref MASSES: HashMap<String, f64> = {
//         let mut m = HashMap::new();
//         m.insert("sun".to_string(), 1.0000000003110439);
//         m.insert("mercury barycenter".to_string(), 0.00000016601);
//         m.insert("venus barycenter".to_string(), 0.0000024478383);
//         m.insert("earth".to_string(), 0.00000300348959632);
//         m.insert("moon".to_string(), 0.00000000007342);
//         m.insert("mars barycenter".to_string(), 0.000000333020);
//         m.insert("jupiter barycenter".to_string(), 0.0009547918932199791); 
//         m.insert("saturn barycenter".to_string(), 0.000285885670706712);
//         m.insert("uranus barycenter".to_string(), 0.00004366249614580214);
//         m.insert("neptune barycenter".to_string(), 0.000051513871954469416);
//         m.insert("pluto barycenter".to_string(), 9.7550000000000000E+02 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000001".to_string(), 6.2628888644409933E+01 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000002".to_string(), 1.3665878145967422E+01 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000003".to_string(), 1.9205707002025889E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000004".to_string(), 1.7288232879171513E+01 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000007".to_string(), 1.1398723232184107E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000010".to_string(), 5.6251476453852289E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000015".to_string(), 2.0230209871098284E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000016".to_string(), 1.5896582441709424E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000031".to_string(), 1.0793714577033560E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000052".to_string(), 2.6830359242821795E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000065".to_string(), 9.3810575639151328E-01 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000087".to_string(), 2.1682320736996910E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000088".to_string(), 1.1898077088121908E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000107".to_string(), 1.4437384031866001E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000433".to_string(), 4.463E-4 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000511".to_string(), 3.8944831481705644E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000704".to_string(), 2.8304096393299849E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m
//     };
// }

lazy_static! {
    pub static ref MASSES: HashMap<String, f64> = {
        let mut m = HashMap::new();
        m.insert("sun".to_string(), 1.3271244004127942E+11 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("mercury barycenter".to_string(), 2.2031868551400003E+04 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("venus barycenter".to_string(), 3.2485859200000000E+05 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("earth".to_string(), 3.9860043550702266E+05 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("moon".to_string(), 4.9028001184575496E+03 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("earth barycenter".to_string(), 3.9860044180000000E+05 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("mars barycenter".to_string(), 4.2828375815756102E+04 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("jupiter barycenter".to_string(), 1.2671276409999998E+08 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("saturn barycenter".to_string(), 3.7940584841799997E+07 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("uranus barycenter".to_string(), 5.7945563999999985E+06 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("neptune barycenter".to_string(), 6.8365271005803989E+06 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("pluto barycenter".to_string(), 9.7550000000000000E+02 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000001".to_string(), 6.2628888644409933E+01 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000002".to_string(), 1.3665878145967422E+01 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000003".to_string(), 1.9205707002025889E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000004".to_string(), 1.7288232879171513E+01 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000007".to_string(), 1.1398723232184107E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000010".to_string(), 5.6251476453852289E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000015".to_string(), 2.0230209871098284E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000016".to_string(), 1.5896582441709424E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000031".to_string(), 1.0793714577033560E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000052".to_string(), 2.6830359242821795E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000065".to_string(), 9.3810575639151328E-01 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000087".to_string(), 2.1682320736996910E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000088".to_string(), 1.1898077088121908E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000107".to_string(), 1.4437384031866001E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000433".to_string(), 4.463E-4 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000511".to_string(), 3.8944831481705644E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000704".to_string(), 2.8304096393299849E+00 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m
    };
}


// BODY2000001_GM = ( 6.2628888644409933E+01 )
//      BODY2000002_GM = ( 1.3665878145967422E+01 )
//      BODY2000003_GM = ( 1.9205707002025889E+00 )
//      BODY2000004_GM = ( 1.7288232879171513E+01 )
//      BODY2000007_GM = ( 1.1398723232184107E+00 )
//      BODY2000010_GM = ( 5.6251476453852289E+00 )
//      BODY2000015_GM = ( 2.0230209871098284E+00 )
//      BODY2000016_GM = ( 1.5896582441709424E+00 )
//      BODY2000031_GM = ( 1.0793714577033560E+00 )
//      BODY2000052_GM = ( 2.6830359242821795E+00 )
//      BODY2000065_GM = ( 9.3810575639151328E-01 )
//      BODY2000087_GM = ( 2.1682320736996910E+00 )
//      BODY2000088_GM = ( 1.1898077088121908E+00 )
//      BODY2000107_GM = ( 1.4437384031866001E+00 )
//      BODY2000433_GM = ( 4.463E-4 )
//      BODY2000511_GM = ( 3.8944831481705644E+00 )
//      BODY2000704_GM = ( 2.8304096393299849E+00 )

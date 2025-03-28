use lazy_static::lazy_static;
use std::collections::HashMap;
use nalgebra::Matrix3;

pub const SPICE_URL: &str = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels";
pub const MPC_URL: &str = "https://www.minorplanetcenter.net/Extended_Files";
pub const MPC_API: &str = "https://data.minorplanetcenter.net/api/get-obs";

pub const KM_TO_AU: f64 = 1.0 / 149_597_870.700;
pub const M_TO_AU: f64 = KM_TO_AU / 1000.0;

pub const SECONDS_PER_DAY: f64 = 86_400.0;

pub const EQUAT_RAD: f64 = 6378137.0;
pub const FLATTEN: f64 = 1.0 / 298.257223563;
pub const O_M_FLATTEN: f64 = 1.0 - FLATTEN;
pub const DEG_TO_RAD: f64 = std::f64::consts::PI / 180.0;

pub const MU_BARY: f64 = 0.00029630927493457475;
// pub const SPEED_OF_LIGHT: f64 = 173.14463268466926; // speed of light in au/day
pub const SPEED_OF_LIGHT: f64 = 173.1446326742403;

pub const GRAVITATIONAL_CONSTANT: f64 = 0.00029591220828559104;
// pub const GRAVITATIONAL_CONSTANT: f64 = 0.00029591220819207774;

pub const ROTATION_J2000: Matrix3<f64> = Matrix3::new(1.0, 0.0, 0.0,
                                                      0.0, 1.0, 0.0,
                                                      0.0, 0.0, 1.0);

pub const ROTATION_FK4: Matrix3<f64> = Matrix3::new(0.999_925_679_495_687_7, 0.011181483239171792, 0.004_859_003_772_314_386,
                                                   -0.01118148322046629, 0.999_937_484_893_313_5, -2.717_029_374_400_203e-5,
                                                   -0.004_859_003_815_359_271, -2.716_259_471_424_704_8e-5, 0.9999881946023742);

pub const ROTATION_GALACTIC: Matrix3<f64> = Matrix3::new(-0.054_875_539_395_742_52, -0.873_437_104_727_596_1, -0.48383499177002515,
                                                         0.494_109_453_627_743_9, -0.44482959429757496, 0.746_982_248_699_891_8,
                                                         -0.867_666_135_683_373_8, -0.19807638961301985, 0.455_983_794_521_419_9);

pub const ROTATION_ECLIPJ2000: Matrix3<f64> = Matrix3::new(1.0, 0.0, 0.0,
                                                           0.0, 0.917_482_062_069_181_8, 0.397_777_155_931_913_7,
                                                           0.0, -0.397_777_155_931_913_7, 0.917_482_062_069_181_8);

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

// lazy_static! {
//     pub static ref MASSES: HashMap<String, f64> = {
//         let mut m = HashMap::new();
//         m.insert("sun".to_string(), 1.327_124_400_412_794_2E11 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("mercury barycenter".to_string(), 2.203_186_855_140_000_3E4 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("venus barycenter".to_string(), 3.248_585_92E5 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("earth".to_string(), 3.986_004_355_070_226_6E5 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("moon".to_string(), 4.902_800_118_457_55E3 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("earth barycenter".to_string(), 3.986_004_418E5 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("mars barycenter".to_string(), 4.282_837_581_575_61E4 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("jupiter barycenter".to_string(), 1.267_127_640_999_999_8E8 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("saturn barycenter".to_string(), 3.794_058_484_18E7 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("uranus barycenter".to_string(), 5.794_556_399_999_998_5E6 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("neptune barycenter".to_string(), 6.836_527_100_580_399E6 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("pluto barycenter".to_string(), 9.755E2 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000001".to_string(), 6.262_888_864_440_993E1 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000002".to_string(), 1.366_587_814_596_742_2E1 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000003".to_string(), 1.920_570_700_202_589 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000004".to_string(), 1.728_823_287_917_151_3E1 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000007".to_string(), 1.139_872_323_218_410_7 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000010".to_string(), 5.625_147_645_385_229 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000015".to_string(), 2.023_020_987_109_828_4 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000016".to_string(), 1.589_658_244_170_942_4 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000031".to_string(), 1.079_371_457_703_356 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000052".to_string(), 2.683_035_924_282_179_5 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000065".to_string(), 9.381_057_563_915_133E-1 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000087".to_string(), 2.168_232_073_699_691 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000088".to_string(), 1.189_807_708_812_190_8 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000107".to_string(), 1.443_738_403_186_6 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000433".to_string(), 4.463E-4 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000511".to_string(), 3.894_483_148_170_564_4 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m.insert("2000704".to_string(), 2.830_409_639_329_985 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
//         m
//     };
// }


lazy_static! {
    pub static ref MASSES: HashMap<String, f64> = {
        let mut m = HashMap::new();
        m.insert("sun".to_string(), 1.327_124_400_412_794_2E11 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("mercury barycenter".to_string(), 2.203_186_855_140_000_3E4 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("venus barycenter".to_string(), 3.248_585_92E5 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("earth".to_string(), 3.986_004_355_070_226_6E5 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("moon".to_string(), 4.9028001184575496E3 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("earth barycenter".to_string(), 4.0350323562548019E+05 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("mars barycenter".to_string(), 4.282_837_581_575_61E4 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("jupiter barycenter".to_string(), 1.267_127_640_999_999_8E8 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("saturn barycenter".to_string(), 3.794_058_484_18E7 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("uranus barycenter".to_string(), 5.794_556_399_999_998_5E6 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("neptune barycenter".to_string(), 6.836_527_100_580_399E6 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("pluto barycenter".to_string(), 9.755E2 * KM3_PER_SECOND2_TO_AU3_PER_DAY2 / GRAVITATIONAL_CONSTANT);
        m.insert("2000001".to_string(), 1.3964518123081070e-13 / GRAVITATIONAL_CONSTANT);
        m.insert("2000002".to_string(), 3.0471146330043200e-14 / GRAVITATIONAL_CONSTANT);
        m.insert("2000003".to_string(), 4.2823439677995011e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000004".to_string(), 3.8548000225257904e-14 / GRAVITATIONAL_CONSTANT);
        m.insert("2000007".to_string(), 2.5416014973471498e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000010".to_string(), 1.2542530761640810e-14 / GRAVITATIONAL_CONSTANT);
        m.insert("2000015".to_string(), 4.5107799051436795e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000016".to_string(), 3.5445002842488978e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000031".to_string(), 2.4067012218937576e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000048".to_string(), 1.9085161956485640e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000052".to_string(), 5.9824315264869841e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000065".to_string(), 2.0917175955133682e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000087".to_string(), 4.8345606546105521e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000088".to_string(), 2.6529436610356353e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000107".to_string(), 3.2191392075878588e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000451".to_string(), 1.2973797046097596e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000511".to_string(), 8.6836253492286545e-15 / GRAVITATIONAL_CONSTANT);
        m.insert("2000704".to_string(), 6.3110343420878887e-15 / GRAVITATIONAL_CONSTANT);
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

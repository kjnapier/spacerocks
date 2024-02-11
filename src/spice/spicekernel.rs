use spice;

pub struct SpiceKernel {
    pub loaded_files: Vec<String>,
}

impl SpiceKernel {

    pub fn new() -> SpiceKernel {
        SpiceKernel { loaded_files: vec![] }
    }

    pub fn load(&mut self, path: &str) {
        spice::furnsh(path);
        self.loaded_files.push(path.to_string());
    }

    pub fn unload(&mut self) {
        spice::kclear();
        self.loaded_files = vec![];
    }

}
use spice;

pub struct SpiceKernel {
    pub loaded_files: Vec<String>,
}

impl SpiceKernel {

    pub fn new() -> SpiceKernel {
        SpiceKernel { loaded_files: vec![] }
    }

    pub fn load(&mut self, path: &str) -> Result<(), String> {
        if self.loaded_files.contains(&path.to_string()) {
            let err = format!("Kernel {} has already been loaded. Skipping.", path);
            return Err(err);
        }
        spice::furnsh(path);
        self.loaded_files.push(path.to_string());
        Ok(())
    }

    pub fn unload(&mut self) {
        spice::kclear();
        self.loaded_files = vec![];
    }

}
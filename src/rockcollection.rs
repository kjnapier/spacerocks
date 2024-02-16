use crate::spacerock::SpaceRock;

struct RockCollection {
    rocks: Vec<SpaceRock>,
}

impl RockCollection {

    pub fn new() -> Self {
        RockCollection {
            rocks: Vec::new(),
        }
    }

    pub fn add(&mut self, rock: SpaceRock) {
        self.rocks.push(rock);
    }

    pub fn remove(&mut self, name: &str) {
        self.rocks.retain(|rock| rock.name != name);
    }

    pub fn get(&self, name: &str) -> Option<&SpaceRock> {
        self.rocks.iter().find(|rock| rock.name == name)
    }

    pub fn get_mut(&mut self, name: &str) -> Option<&mut SpaceRock> {
        self.rocks.iter_mut().find(|rock| rock.name == name)
    }
    
}

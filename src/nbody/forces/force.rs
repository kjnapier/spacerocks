use crate::spacerock::SpaceRock;
use nalgebra::Vector3;

pub trait Force: Send + Sync + ForceClone {
    fn apply(&self, entities: &mut Vec<SpaceRock>);
}


trait ForceClone {
    fn clone_box(&self) -> Box<dyn Force + Send + Sync>;
}

impl<T> ForceClone for T
where
    T: 'static + Force + Clone,
{
    fn clone_box(&self) -> Box<dyn Force + Send + Sync> {
        Box::new(self.clone())
    }
}

// We can now implement Clone manually by forwarding to clone_box.
impl Clone for Box<dyn Force + Send + Sync>{
    fn clone(&self) -> Box<dyn Force + Send + Sync> {
        self.clone_box()
    }
}

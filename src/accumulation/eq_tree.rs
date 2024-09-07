use ark_ff::PrimeField;


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EqTree<F: PrimeField> {
    // vector of length 2 * 2^depth - 1
    pub nodes: Vec<F>,
    pub depth: usize,
}

impl<F: PrimeField> EqTree<F> {

    /// generates the eq tree given a vector of length depth
    pub fn new(x: &[F]) -> Self {
        let depth = x.len();
        let mut nodes = vec![F::ZERO; 2 * (1 << depth) - 1];

        // The root node starts with value 1.
        nodes[0] = F::ONE;

        // Build the tree
        for i in 0..depth {
            for j in 0..(1 << i) {
                let node_idx = (1 << i) + j - 1;
                let val = nodes[node_idx];
                nodes[2 * node_idx + 1] = val * (F::ONE - x[i]);
                nodes[2 * node_idx + 2] = val * x[i];
            }
        }

        Self { nodes, depth }
    }

    /// generates the error values for a tree given a vector
    pub fn difference(&self, x: &[F]) -> Self {
        assert_eq!(x.len(), self.depth, "inconsistent depth");
        let depth = x.len();
        let mut nodes = vec![F::ZERO; 2 * (1 << depth) - 1];

        // The root node starts with value 0.
        nodes[0] = F::ZERO;

        // Build the tree
        for i in 0..depth {
            for j in 0..(1 << i) {
                let node_idx = (1 << i) + j - 1;
                println!("{}", node_idx);
                let val = self.nodes[node_idx];
                nodes[2 * node_idx + 1] = self.nodes[2 * node_idx + 1] - val * (F::ONE - x[i]);
                nodes[2 * node_idx + 2] = self.nodes[2 * node_idx + 2] - val * x[i];
            }
        }

        Self { nodes, depth }
    }

    /// checks all values in a tree are zero
    pub fn is_zero(&self) -> () {
        for i in &self.nodes {
            assert!(i.is_zero());
        }
    }

    /// returns leaves of the tree starting from (1-x_1)...(1-x_n) to x_1...x_n
    pub fn get_leaves(&self) -> &[F] {
        &self.nodes[(1 << (self.depth)) - 1..]
    }

    /// prints the different layers of the tree one by one
    pub fn print_tree(&self) {
        let mut level_start = 0;
        for i in 0..self.depth {
            let num_nodes_at_level = 1 << i;
            let level_end = level_start + num_nodes_at_level;
            let level_nodes = &self.nodes[level_start..level_end];

            // Print the depth and the nodes at this level
            print!("Depth {}: ", i);
            for node in level_nodes {
                print!("{:?} ", node);
            }
            println!();

            level_start = level_end;
        }

        // Printing leaves (depth `depth`)
        print!("Depth {}: ", self.depth);
        for leaf in &self.nodes[level_start..] {
            print!("{:?} ", leaf);
        }
        println!();
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::{AdditiveGroup, Field};

    use crate::constant_for_curves::ScalarField;

    use super::*;

    type F = ScalarField;

    #[test]
    fn test_tree() {
        let x = vec![F::ONE, F::ZERO];
        let tree = EqTree::new(x.as_slice());
        tree.print_tree();
        let dif = tree.difference(x.as_slice());
        dif.is_zero();
    }
}

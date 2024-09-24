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
        // reversing the array
        let x = {
            let mut temp = x.to_vec();
            temp.reverse();
            temp
        };

        // define depth the length of the array
        let depth = x.len();

        // initialise all nodes as zero
        let mut nodes = vec![F::ZERO; 2 * (1 << depth) - 1];

        // The root node starts with value 1.
        nodes[0] = F::ONE;

        // Build the tree
        for i in 0..depth {
            for j in 0..(1 << i) {
                let node_idx = (1 << i) + j - 1;
                let val = nodes[node_idx];
                let left_index = (2 * (1 << i) - 1) + j;
                let right_index = left_index + (1 << i);
                nodes[left_index] = val * (F::ONE - x[i]);
                nodes[right_index] = val * x[i];
            }
        }

        Self { nodes, depth }
    }

    /// generates the error values for a tree given a vector
    pub fn difference(&self, x: &[F]) -> Self {
        assert_eq!(x.len(), self.depth, "inconsistent depth");

        // reversing the array
        let x = {
            let mut temp = x.to_vec();
            temp.reverse();
            temp
        };

        let depth = x.len();
        let mut nodes = vec![F::ZERO; 2 * (1 << depth) - 1];

        // The root node starts with value 0.
        nodes[0] = F::ZERO;

        // Build the tree
        for i in 0..depth {
            for j in 0..(1 << i) {
                let node_idx = (1 << i) + j - 1;
                let val = self.nodes[node_idx];
                let left_index = (2 * (1 << i) - 1) + j;
                let right_index = left_index + (1 << i);
                nodes[left_index] = self.nodes[left_index] - val * (F::ONE - x[i]);
                nodes[right_index] = self.nodes[right_index] - val * x[i];
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
    use ark_std::UniformRand;
    use rand::thread_rng;

    use crate::constant_for_curves::ScalarField;
    use crate::polynomial::eq_poly::EqPolynomial;

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

    #[test]
    fn test_tree2() {
        let x = vec![
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng()),
            F::rand(&mut thread_rng()),
        ];
        let tree = EqTree::new(x.as_slice());
        let results: Vec<F> = EqPolynomial::new(x).evals();

        assert_eq!(tree.get_leaves().to_vec(), results);
    }
}
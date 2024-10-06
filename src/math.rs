// borrowed from Nexus

pub trait Math {
    fn square_root(self) -> usize;
    fn pow2(self) -> usize;
    fn get_bits_non_canonical_order(self, num_bits: usize) -> Vec<bool>;
    fn get_bits_canonical_order(self, num_bits: usize) -> Vec<bool>;
    fn log_2(self) -> usize;
}

impl Math for usize {
    #[inline]
    fn square_root(self) -> usize {
        (self as f64).sqrt() as usize
    }

    #[inline]
    fn pow2(self) -> usize {
        let base: usize = 2;
        base.pow(self as u32)
    }

    /// Returns the num_bits from n in a non-canonical order
    fn get_bits_non_canonical_order(self, num_bits: usize) -> Vec<bool> {
        (0..num_bits)
            .map(|shift_amount| (self & (1 << (num_bits - shift_amount - 1))) > 0).rev()
            .collect::<Vec<bool>>()
    }

    /// Returns the num_bits from n in a canonical order
    fn get_bits_canonical_order(self, num_bits: usize) -> Vec<bool> {
        (0..num_bits)
            .map(|shift_amount| (self & (1 << (num_bits - shift_amount - 1))) > 0)
            .collect::<Vec<bool>>()
    }

    fn log_2(self) -> usize {
        assert_ne!(self, 0);

        if self.is_power_of_two() {
            (1usize.leading_zeros() - self.leading_zeros()) as usize
        } else {
            (0usize.leading_zeros() - self.leading_zeros()) as usize
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::math::Math;

    #[test]
    fn test() {
        println!("{:?}", 11.get_bits_non_canonical_order(8));
        assert_eq!(2000f64.log2().floor() as usize, 10);
        assert_eq!(10f64.log2().floor() as usize, 3);
        assert_eq!(1024f64.log2().floor() as usize, 10);
        assert_eq!(0f64.log2().floor() as usize, 0);
    }
}

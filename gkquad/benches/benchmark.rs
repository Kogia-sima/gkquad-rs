extern crate gkquad;

mod double;
mod single;

use double::double;
use single::single;

use smbench::*;

smbench_trace_memory!();
smbench_main!(single, double);

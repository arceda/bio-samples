
<?php

echo read_fasta("uploads/POLSPEVP1.fasta");

function read_fasta($target_file){     

     $f = file_get_contents($target_file);
     $array = explode("\n",$f);
     unset($array[0]);

     //var_dump($array);
     $seq = implode('', $array);

     return trim($seq);
}

?>
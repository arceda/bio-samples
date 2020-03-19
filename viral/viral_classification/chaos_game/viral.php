<?php
header("Access-Control-Allow-Origin: *");

$model =  $_POST["model"];

$fileName = $_FILES['file']['name'];
$fileType = $_FILES['file']['type'];
$fileError = $_FILES['file']['error'];
$fileContent = file_get_contents($_FILES['file']['tmp_name']);

if($fileError == UPLOAD_ERR_OK){

    $array = explode('.', basename($fileName));
	$extension = end($array);
    $name = $fileName. "_".uniqid(rand(), true) . microtime(true) . "." . $extension;
	$target_dir = "/var/www/html/viral/uploads/";
	$target_file = $target_dir . $name;


	if (move_uploaded_file($_FILES["file"]["tmp_name"], $target_file)) {
          $output = process($target_file, $model);
          $seq = read_fasta($target_file);

          $output = str_replace('[', '', $output);
          $output = str_replace(']', '', $output);
          $output = str_replace('"', '', $output);
          $output = str_replace("'", '', $output);

		$response = array('result' => 1, 'message' => $output, 'seq' => $seq);
		echo json_encode($response);

		return;

	} else {
		$response = array('result' => 0, 'message' => 'Problem uploading', 'seq' => '');
		echo json_encode($response);
	}
   

}else{
   switch($fileError){
     case UPLOAD_ERR_INI_SIZE:   
          $message = 'Error al intentar subir un archivo que excede el tamaño permitido.';
          break;
     case UPLOAD_ERR_FORM_SIZE:  
          $message = 'Error al intentar subir un archivo que excede el tamaño permitido.';
          break;
     case UPLOAD_ERR_PARTIAL:    
          $message = 'Error: no terminó la acción de subir el archivo.';
          break;
     case UPLOAD_ERR_NO_FILE:    
          $message = 'Error: ningún archivo fue subido.';
          break;
     case UPLOAD_ERR_NO_TMP_DIR: 
          $message = 'Error: servidor no configurado para carga de archivos.';
          break;
     case UPLOAD_ERR_CANT_WRITE: 
          $message= 'Error: posible falla al grabar el archivo.';
          break;
     case  UPLOAD_ERR_EXTENSION: 
          $message = 'Error: carga de archivo no completada.';
          break;
     default: $message = 'Error: carga de archivo no completada.';
              break;
    }
      echo json_encode(array(
               'result' => 0,
               'message' => $message,
               'seq' => ''
            ));
}



function process($target_file, $model){
    $exe = "/home/vicente/projects/BIOINFORMATICS/bio-samples/viral/viral_classification/test.py '$target_file' $model";

    //echo "python '$exe' '$first_name' '$path_image'";

    exec("python3 $exe", $output);

    return $output;
}

function read_fasta($target_file){     

     $f = file_get_contents($target_file);
     $array = explode("\n",$f);
     unset($array[0]);

     //var_dump($array);
     $seq = implode('', $array);

     return trim($seq);
}

?>
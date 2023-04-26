#!/bin/bash


# Binary data of the personalized picture
image_data=$(xxd -p -c 99999 Eagloo_frontend_functions/eagloo_logo_drawing.png | tr -d '\n')

# Write the personalized picture to a file
cat << EOF | xxd -r -p > personalized_image.png
$image_data
EOF

# Display the personalized picture
feh personalized_image.png

source ./dependencies/bin/activate

python ./Eagloo_frontend.py

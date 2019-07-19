const path = require("path"),
  MiniCssExtractPlugin = require('mini-css-extract-plugin');

module.exports = {
  entry: path.resolve("viz", "app.jsx"),
  output: {
    path: path.resolve(__dirname, "convex_qsr", "static"),
    filename: "main.js"
  },
  plugins: [
    new MiniCssExtractPlugin({
      filename: 'style.css',
    })
  ],
  module: {
    rules: [
      {
        test: /\.jsx$/,
        exclude: /node_modules/,
        loader: "babel-loader",
        query: {
          presets: ["@babel/react"]
        }
      },
      {
        test: /\.(sa|sc|c)ss$/,
        use: [
          {
            loader: MiniCssExtractPlugin.loader,
          },
          'css-loader',
          'postcss-loader',
          'sass-loader',
        ],
      }
    ]
  },
  mode: "development",
  devtool: "inline-source-map"
};

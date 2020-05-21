const path = require("path"),
  HtmlWebpackPlugin = require("html-webpack-plugin"),
  UglifyJsPlugin = require("uglifyjs-webpack-plugin"),
  MiniCssExtractPlugin = require("mini-css-extract-plugin");


const static_dir = path.resolve(__dirname, "superseal", "viz"),
  port = 8123;

module.exports = {
  entry: path.resolve("viz", "app.jsx"),
  output: {
    path: static_dir,
    filename: "main.js"
  },
  plugins: [
    new HtmlWebpackPlugin({
      title: "SuperSEAL"
    }),
    new MiniCssExtractPlugin({
      filename: "style.css",
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
          "css-loader",
          "postcss-loader",
          "sass-loader",
        ],
      }
    ]
  },
  mode: "development",
  devtool: "inline-source-map",
  devServer: {
    port: port,
    historyApiFallback: true,
    contentBase: static_dir
  },
  optimization: {
    minimizer: [new UglifyJsPlugin()]
  }
};
